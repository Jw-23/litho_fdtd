# litho_fdtd

`litho_fdtd` 是一个 C++23 三维时域有限差分（FDTD）光刻掩模频谱求解器。它复用
`litho_geometry` 的几何类型，补齐 Yee 网格 Maxwell 求解器，输出可直接用于成像模型或
光源掩模优化的复振幅衍射级次。

公共命名空间为 `litho::fdtd`。接口风格与 `litho_geometry` 一致：值语义、自由函数、
`std::expected` 错误返回、不显式 `throw`，源码注释使用英文。

## 功能

- 三维笛卡尔 Yee 网格，内部二阶差分。
- 任意非周期方向使用卷积 PML（CPML）吸收边界。
- 横向周期与 Bloch 周期边界，支持掩模周期和斜入射。
- 在均匀层中初始化的脉冲平面波光源；对 s/p 偏振都使用一致的 E/H 初值和包络导数修正。
- 由参考波长复折射率构造各向同性材料：正实部使用损耗介质模型，负实部使用单 Drude 极点。
- 平面 DFT 监视器和掩模频谱后处理：复振幅级次、基于真实场 Poynting 矢量的效率、
  归一化频谱。
- 通过 `litho_geometry::common::extrusion` 和 `box3` 描述三维掩模；逐 z 切片体素化。
- SoA 存储、融合的 H/E 调度 phase，以及来自同级 `parallel` 项目的分块调度器；场更新最多使用
  4 个性能核。

## 构建

依赖：CMake >= 3.25、Ninja（推荐）、支持 `std::expected` 的 C++23 编译器，Boost 头文件，以及同级的 `parallel` 调度器库。若尚未定义 `parallel` 目标，
`litho_fdtd` 会自动添加 `../parallel`。

```sh
cmake --workflow --preset release
./build/release/tests/litho_fdtd_tests
./build/release/examples/litho_fdtd_mask_example
```

也可以分步执行：

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

如果 `litho_geometry::litho_geometry` 尚未定义，`litho_fdtd` 会自动添加同级的
`litho_geometry` 子项目。将两个项目一起纳入上层工程时，可以直接：

```cmake
add_subdirectory(path/to/litho_fdtd)
target_link_libraries(my_app PRIVATE litho_fdtd::litho_fdtd)
```

构建吞吐基准：

```sh
cmake --preset release -DLITHO_FDTD_BUILD_BENCHMARKS=ON
cmake --build --preset release
./build/release/benchmarks/litho_fdtd_benchmark [workers]
```

基准会输出 120×120×240 真空网格的百万单元更新/秒；可选参数对应 `run_options::workers`。

## 快速开始：193 nm 线/空掩模

```cpp
#include <litho_fdtd/fdtd.hpp>

using namespace litho::fdtd;
namespace common = litho::geometry::common;

mask_problem p;
p.grid.cells = {80, 4, 120};
p.grid.spacing = {5e-9, 5e-9, 5e-9};

// add() 返回 std::expected<size_t>；示例假定成功。
(void)p.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt});
(void)p.materials.add(*from_refractive_index("fused silica", {1.56, 0.0}, 193e-9));
(void)p.materials.add(*from_refractive_index("MoSi", {2.5, 0.6}, 193e-9));

p.mask.background = 0;
p.mask.layers.push_back(planar_layer{0.0, 300e-9, 1}); // 石英衬底

auto line = common::make_extrusion(
    *common::rectangle(common::box<double, 2>{{0.0, 0.0}, {200e-9, 20e-9}}),
    300e-9, 360e-9);
p.mask.solids.push_back(patterned_solid{*line, 2});

p.source.z = 120e-9;         // 初始脉冲中心，位于石英内部
p.source.wavelength = 193e-9;
p.source.polarization = 0;   // 0 为 s（TE），pi/2 为 p（TM）
p.source.pulse_fwhm = 0.5;   // 光周期
p.run.duration = 20;         // 光周期
p.run.workers = 0;           // 0 = 4 个性能核（最多 4）
p.run.chunk_rows = 128;      // 调度器行粒度

p.monitor_z = 450e-9;        // 空气中，吸收层上方
p.spectrum.min_order_x = -6;
p.spectrum.max_order_x = 6;
p.spectrum.background_index = 1.0;

auto result = compute_mask_spectrum(p);
if (!result) { /* result.error(), message(result.error()) */ }

for (const auto& order : result->spectrum.orders) {
    // order.field 为 {Ex, Ey, Ez}；有参考时 order.efficiency 有值。
    // result->normalized 为除以参考 (0,0) 级次后的复振幅。
}
```

可运行示例为 [`examples/mask_spectrum.cpp`](examples/mask_spectrum.cpp)，会打印各级次并写出
`litho_fdtd_mask_spectrum.csv`。

## 三维几何

三维定义位于 `litho_geometry/common/solid.hpp`，并通过 `<litho_geometry/common.hpp>` 导出：

- `common::extrusion<T>`：xy 多边形沿 z 从 `z_min` 拉伸到 `z_max`，`top_scale` 可描述
  线性收窄的侧壁。
- `common::box3<T>` / `box3d`：轴对齐三维长方体。
- `common::solid3<T>` / `solid3d`：`std::variant<extrusion<T>, box3<T>>`。
- `common::section_at(solid, z)`：返回精确二维截面，供逐切片体素化使用。
- `common::bounds(solid)`：三维包围盒。

求解器逐 z 切片光栅化，不会在几何层构造完整体素体。实体覆盖平面层，后加入的实体覆盖先加入的
实体。

## 接口概览

| 头文件 | 内容 |
| --- | --- |
| `<litho_fdtd/grid.hpp>` | `grid`、`axis`、`component`、Yee 偏移 |
| `<litho_fdtd/material.hpp>` | `material`、`material_library`、`from_refractive_index`、`from_permittivity` |
| `<litho_fdtd/structure.hpp>` | `planar_layer`、`patterned_solid`、`structure`、`voxelize` |
| `<litho_fdtd/boundary.hpp>` | `boundaries`、`boundary_kind`、`cpml_options` |
| `<litho_fdtd/source.hpp>` | `plane_wave_source` |
| `<litho_fdtd/monitor.hpp>` | `monitor`、`monitor_result`、平面 DFT 累积 |
| `<litho_fdtd/spectrum.hpp>` | `spectrum_options`、`diffraction_order`、`compute_spectrum`、`normalize_spectrum` |
| `<litho_fdtd/solver.hpp>` | `problem`、`run_options`、`run_result`、`solve` |
| `<litho_fdtd/mask.hpp>` | `mask_problem`、`mask_spectrum`、`compute_mask_spectrum` |
| `<litho_fdtd/fdtd.hpp>` | 汇总头文件 |

所有可能失败的入口返回 `litho::fdtd::result<T>`（即 `std::expected<T, errc>`）。实现只在
公共边界捕获分配异常并转换为 `errc::allocation_failure`，调用方不需要处理异常。

## 单位与约定

- 内部使用 SI 单位：`E` 为 V/m，`H` 为 A/m，长度为 m，时间为 s。
- 时间因子为 `exp(-i omega t)`。复折射率输入为 `n = n' + i n''`，无源材料要求 `n'' >= 0`。
- `from_refractive_index` 根据 `Re(n^2)` 选择：
  - 正实部：非色散损耗介质，`epsilon_r = Re(n^2)`，
    `sigma = omega * eps0 * Im(n^2)`；
  - 负实部：单 Drude 极点，在参考波长精确匹配 `n^2`，从而保持时域更新无源、稳定，
    可表示金属和负介电常数吸收层。
- 角度均为弧度。`theta` 为相对 `+z` 的入射角，`phi` 为 xy 平面方位角。
  `polarization = 0` 为 s（TE），`pi/2` 为 p（TM）。`from_below` 选择 +z 或 -z 传播。
- `plane_wave_source::z` 是初始脉冲中心，不是电流片。脉冲是单向波包，不会产生反向传播的
  源波。

## 数值方法

### Yee 网格

六个分量共用同一数组形状，并按单元下棱/下面交错：

| 分量 | 位置 |
| --- | --- |
| `Ex` | `(i+1/2, j, k)` |
| `Ey` | `(i, j+1/2, k)` |
| `Ez` | `(i, j, k+1/2)` |
| `Hx` | `(i, j+1/2, k+1/2)` |
| `Hy` | `(i+1/2, j, k+1/2)` |
| `Hz` | `(i+1/2, j+1/2, k)` |

E 更新对 H 使用后向差分，H 更新对 E 使用前向差分。非周期外边界处最外层使用单侧差分，
CPML 在波到达边界前吸收。

### 时间步长

`dt = courant / (c_max * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2))`，其中 `c_max` 由材料
`epsilon_inf * mu_r` 的最小值估计最大相速。默认 Courant 数为 0.95。含 Drude 材料时还会检查
`omega_p * dt <= 1`。

### CPML

每个非周期方向都有多项式渐变 CPML，并分别提供整数/半整数节点系数表，以匹配交错场位置。
默认 `layers = 10`、`reflection = 1e-6`、`order = 3`、`kappa_max = 5`、
`alpha_max = 0.05 * sigma_max`。光学掩模建议 PML 厚度达到半波长以上；5–10 nm 间距下
通常使用 20–40 层。PML 对几何透明，不需要在结构中额外处理。

### 周期与 Bloch 边界

`boundary_kind::periodic` 的方向在更新模板中环绕。斜入射时求解器自动切换到复解析场，并在
环绕处施加 Bloch 相位 `exp(i k_parallel L)`。正入射或横向为 CPML 时使用实场，节省一半内存
和计算量。

### 光源归一化

脉冲在 `t = 0`（E）和 `t = -dt/2`（H）初始化为高斯平面波包。H 初值包含包络导数项，使初始
条件对脉冲平面波满足 Maxwell 方程，而不仅对单频波成立。监视器 DFT 累加会除以解析高斯包络
积分，因此单位振幅平面波给出单位复频谱系数。`monitor::time_window` 可用于排除 CPML 或腔体
后到回波；此时归一化不再是完整脉冲积分，比较完整与参考运行时请使用相同窗口。

### 掩模频谱

z 法向监视器的级次约定为：

```text
kx = kx_incident + 2*pi*m/period_x
ky = ky_incident + 2*pi*n/period_y
kz = sqrt((n_bg*k0)^2 - kx^2 - ky^2)
```

后处理读取平面上的六个 DFT 场分量，使用 E、H 的真实 Fourier 系数计算级次功率
（`S_z = 0.5 Re(Ex Hy* - Ey Hx*)`），并补偿 E 平面与交错 H 采样之间的半单元空间相位。
倏逝级次报告零传播功率。反射监视器需设置 `propagation_sign = -1`；默认 `+1` 为透射。

`compute_spectrum` 返回复振幅级次；设置 `incident_power` 后同时返回效率。
`normalize_spectrum(full, reference)` 用参考 (0,0) 级次的主分量归一化所有级次并保留相对相位。
`compute_mask_spectrum` 在 `use_reference = true` 时自动执行完整结构与无图形参考叠层两次
仿真，返回原始级次、归一化级次、入射功率和总透射效率。

## 精度建议

- 最高折射率材料中的最短波长至少 20 个网格；相位敏感的掩模频谱建议 30–40 个网格。
- CPML 厚度按波长计要足够。光学波段十层 PML 可能反射几个百分点；λ/20 间距下 20–40 层更稳妥。
- 将初始脉冲放在均匀、无损耗层中，并使脉冲包络远离掩模和 PML。求解器会检查光源周围
  ±4σ_z 的切片，拒绝跨越材料边界或进入 CPML 的光源。
- 材料界面在交错场位置光栅化，平面界面存在等效半单元偏移。关键界面可对齐到期望的半单元
  平面，或细化网格。当前没有亚像素界面平滑。
- 效率应使用参考运行。高层接口在 `use_reference = true` 时会这样做；参考运行会移除全部
  `solids`，因此衬底等无图形层应放在 `mask.layers`。
- 材料模型在参考波长精确。定量频谱请使用围绕该波长的窄带脉冲；宽带问题请按波长分别运行，
  或拟合多极点模型。

## 验证

`tests/` 包含：

- 三维拉伸/实体切片、包围盒与体素化。
- 复折射率转换、Drude 往返、主动材料拒绝。
- 真空平面波振幅、波阻抗与能量。
- 正入射 Fresnel 反射与 `(1-n)/(1+n)` 对比。
- 斜入射 s/p 偏振 Bloch 周期平面波，包括横向 `E · k`。
- 无损介质光栅能量守恒（`T + R ≈ 1`，使用真实场级次功率）。
- 高层掩模频谱接口与低层光栅结果的一致性。
- Drude 吸收层运行，输出有界且有限。

开发机上 Release 测试套件约 4–6 秒通过；Debug 因更新核未优化会明显更慢。

## 性能与内存

求解器存储六个场数组、最多十二个 CPML 辅助数组（仅 CPML 方向需要）以及三个电场和三个磁场
的对角材料系数数组。复场（Bloch）运行时场和 CPML 数组使用 `std::complex<double>`。

更新核通过同级 `parallel` 工作窃取调度器的 `parallel::for_each_chunk` 分发。调度器把
z×y 行切成适合缓存的块；调用线程参与计算，同一线程的空闲线程组会在多次更新之间复用。
`run_options::workers = 0` 选择 4 个性能核，`1` 为串行，显式更大的值也会被限制为 4，
避免 macOS 能效核进入场更新。`run_options::chunk_rows` 指定调度器以行为单位的粒度
（默认 128；`0` 使用框架默认值）。监视器累积和材料设置仍为单线程，因为它们相对场更新很小。
不依赖 OpenMP 运行时或构建选项。

H 和 E 更新融合为每个时间步两个调度 phase：三个 H 分量共用一个并行区域，随后三个 E 分量
共用一个并行区域。CPML 辅助变量只在 PML 层内读写；在 PML 外直接使用普通 Yee 差分。

开发机上，基准在 Release、120×120×240 网格、300 步、193 nm 真空下测得：1 worker 约
130 Mcell-updates/s，4 worker 约 344，加速比约 2.6 倍。大掩模通常受共享内存带宽限制；
网格控制在脉冲和 CPML 所需的合理范围即可。

## 限制

- 频谱提取仅支持 z 法向监视器；仅记录场的监视器可使用 x 或 y 法向。
- 只有一个 Drude 极点，没有多极点 Lorentz/Drude、各向异性或磁损耗。
- 没有 MPI 域分解，只有 `parallel` 的共享内存调度器。
- 没有亚像素界面平滑，也没有共形 FDTD。
- 高层参考运行会移除全部实体；图形化衬底应拆成平面层加独立三维实体。

## 相关文档

- [English documentation](README.md)
- [`litho_geometry` 几何库](../litho_geometry/README.md)

## 配置检查与监视器定位

`validate(problem)` 可在体素化、场分配之前检查网格、材料引用、边界、光源参数、运行参数和
监视器。依赖体素化结果的光源均匀性和材料相关稳定性仍由 `solve` 检查。
显式设置 `run.steps` 时，`duration` 不参与校验或计算。

```cpp
auto detector = monitor_at(p.grid, axis::z, 450e-9);
if (!detector) { /* message(detector.error()) */ }
// 成功后可设置 name、frequencies、time_window、spectrum，并加入 problem.monitors。
// compute_mask_spectrum 内部也使用相同的坐标转换。
```

`monitor_at` 支持 x/y/z 法向，按 SI 物理坐标定位；允许的浮点偏差为 1e-6 个网格，
不在网格平面上的坐标会返回错误。`component_values` 对无效频率、分量和不完整缓冲区返回空 span。
`run_result::workers` 现在报告实际选择的线程数，自动模式的小网格会报告 1。

## 日志与性能分析

构建会复用已有的 `litho::trace` 目标，否则从 `../litho_tracing` 构建采集库；不会构建其分析工具、
示例或测试。可用 `-DLITHO_FDTD_TRACING_SOURCE_DIR=/path/to/litho_tracing` 指定其他源码路径，
也可在上层先 `find_package(litho_trace CONFIG REQUIRED)`。

```cpp
litho::Options options;
options.directory = "runs/fdtd-001"; // 每次运行使用新目录
options.build_id = "fdtd-release";
options.source_revision = "your-git-revision";
auto opened = litho::Session::open(options);
if (!opened) { /* opened.error().describe() */ }
auto& trace = opened.value();

p.run.trace.session = &trace; // problem 和 mask_problem 都通过 run 配置
p.run.trace.context = trace.stream("mask-job", {}, LT_SOURCE);
p.run.trace.step_interval = 1;   // 每步采样；0 关闭 kernel span
p.run.trace.detail_steps = 1;    // 第 0 步额外记录 per-chunk worker span
{
    auto result = compute_mask_spectrum(p);
    if (!result) { /* message(result.error()) */ }
}
trace.close();
auto stats = trace.stats(); // 检查 dropped_events 和 error
```

Session 由调用方管理，必须覆盖求解调用的整个生命周期。默认空指针不记录；求解器不创建目录、
不输出控制台日志。插桩按尽力方式工作，同步采集异常不会改变求解结果；正常采集后的丢弃与写盘错误
可通过 Session 统计检查。不要在另一个线程中销毁仍被求解器使用的 Session。

- `fdtd`：`mask_spectrum`、`reference` / `patterned`、`solve`、`voxelize`、`coefficients`、
  `run`、`initialize`、`time_loop`、`postprocess` 阶段，以及配置、准备完成、完成和失败事件。
- `fdtd.step`：每组时间步的 span 和进度，包含起始步、已完成步和总步数。
- `fdtd.kernel`：每组**第一步**记录融合的 `update_h` / `update_e` 和 monitor span。每个
  phase span 带一个 `summary` instant：
  - `wall_ns`：调用线程看到的 phase 墙钟时间；
  - `main_work_ns`：调用线程自己执行 chunk 的时间；
  - `wait_ns`：`wall_ns - main_work_ns`；
  - `work_ns`：所有 worker chunk 时长之和（仅 detail 步）；
  - `chunks`：实际 chunk 数；
  - `utilization`：`work_ns / (wall_ns * workers)`（仅 detail 步）；
  - `main_fraction`：`main_work_ns / wall_ns`。
- `fdtd.chunk`：`detail_steps` 指定前 N 个采样步的 worker chunk span，用于在时间线中查看各
  worker 的活跃区间和 barrier 尾部。`utilization` 接近 1 但加速比低，说明单 worker 吞吐下降
  （内存/缓存竞争）；`utilization` 低说明 worker 空闲或负载不均。可通过 `Options::categories`
  过滤类别。

事件包含源码位置与显式父子 span；参考运行和完整结构运行可以分别查看。
`elapsed_seconds` 仍只计时间循环（启用采集时包含循环内插桩开销）；外层 `solve` span 包含准备和
后处理，基准还输出整个 solve 的墙钟时间。工作线程内不逐单元记录事件。

```sh
./build/release/benchmarks/litho_fdtd_benchmark 4 runs/fdtd-w4
./build/release/benchmarks/litho_fdtd_benchmark 4 monitors runs/fdtd-benchmark-monitors
../litho_tracing/build/litho_analyze runs/fdtd-w4 summary
../litho_tracing/build/litho_analyze runs/fdtd-w4 serve web/dist 8090
../litho_tracing/build/litho_analyze runs/fdtd-w4 export runs/fdtd-w4/perfetto.json
./build/release/benchmarks/litho_fdtd_spectrum_benchmark
```

分析工具需在 `litho_tracing` 中单独构建。`monitors` 基准添加三个 z 平面，每个平面记录四个频率。
频谱基准使用 120×120 平面、169 个级次，重复五次并打印校验值。
频谱提取将复指数按 x/y 分离并提出交错偏移，把每级次的复指数求值从 `5*nx*ny` 降到
`nx+ny+2`；监视器按频率/分量连续遍历，不再在每一步分配相位数组。

本机 AppleClang 21、Release、4 worker 的前后对照（相同任务串行运行，未开启采集）：

| 场景 | 优化前 | 优化后 | 变化 |
| --- | ---: | ---: | --- |
| 120×120 平面、169 级次，五次平均 | 81.023 ms | 13.703 ms | 约 5.9 倍速度 |
| 120×120×240、300 步、三个四频监视器，三次时间循环中位数 | 4.553 s | 4.254 s | 耗时减少 6.6% |

在相同的 120×120×240、300 步真空网格、无监视器条件下（`step_interval = 1`、
`detail_steps = 1`），对优化前后做 4 worker 对比：

| 配置（4 workers） | 时间循环 | 第 0 步 phase 总工作量 | worker utilization | worker 空闲 |
| --- | ---: | ---: | ---: | ---: |
| 优化前：6 phase、全网格 CPML、64 行粒度 | 4.63 s | 57.1 ms | 0.805 | 19.5% |
| H/E 融合、边界 CPML、64 行粒度 | 3.27 s | 40.9 ms | 0.924 | 7.6% |
| 同上，128 行粒度（默认） | 3.28 s | 44.0 ms | 0.956 | 4.4% |
| 同上，256 行粒度 | 3.29 s | 43.5 ms | 0.969 | 3.1% |
| 同上，512 行粒度 | 3.36 s | 42.5 ms | 0.976 | 2.4% |

单线程时间循环从优化前 12.11 s 降到约 8.0–8.3 s。4 worker 第 0 步加速比为
`4 * 0.956 / 1.56 = 2.45`，与实测约 2.5 倍一致：融合和边界 CPML 降低了绝对工作量，但并发
worker 仍因内存/缓存带宽和单核吞吐下降使总工作量增大到 1.56 倍。增大粒度会降低 worker 空闲
（19.5% → 4.4% → 2.4%），但尾部变粗后墙钟时间不再改善；128 行是默认平衡点。主线程自身只
等待几个百分点，之前的等待主要发生在 worker 的 phase 边界。三次采集写入 3,622–5,430 条
事件、零丢弃。数据与硬件相关。
后续可优先评估批量任务中的参考频谱复用，以及网格/材料系数的预编译与复用；当前高层接口仍按每次
请求重新体素化并执行参考仿真，尚未加入缓存生命周期接口。
