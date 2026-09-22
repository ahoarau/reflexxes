# reflexxes

Reflexxes Type IV online trajectory generation library, copied from
[intrinsic-core](https://github.com/intrinsic-ai) (`intrinsic_control/intrinsic/icon/reflexxes`)
with a standalone CMake build. Apache-2.0 (see `LICENSE`).

Sources keep their original include paths (`intrinsic/icon/reflexxes/...`, `intrinsic/util/...`)
under `src/`. The only external dependency is Abseil.

```sh
pixi run build                    # configure + build
pixi run build --target reflexxes # build a single target
pixi run test
pixi run example
```

Use from another CMake project:

```cmake
find_package(reflexxes CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE reflexxes::reflexxes)
```
