# SQPhysics2D

SQPhysics2D is an experimental C++20 2D physics engine that advances connected
groups of bodies, called physics islands, on their own schedules. Its
timestamped command queue and in-memory snapshots let you restore a simulation
and check its replay one tick at a time.

```text
Timestamped commands
        |
        v
Schedule physics islands
        |
        v
Detect collisions and solve contacts
        |
        v
Record a state hash at each complete tick
```

## Replay demo

The [example](examples/replay.cpp) creates a rigid body with a circle shape,
captures its world, and schedules an impulse two ticks ahead. It then restores
the snapshot and submits the same input. If the runs diverge, it reports the
first differing tick.

From a clone, use CMake 3.24+ and a C++20 compiler. These commands use Windows
with Visual Studio and MSVC.

```powershell
git submodule update --init
cmake -S . -B build -A x64
cmake --build build --config Release --parallel
.\build\bin\Release\sqphysics2d_replay.exe
ctest --test-dir build -C Release --output-on-failure
```

```text
Replay matched 8 ticks.
```

State hashes are comparable within the same build and process. When replaying
a world that uses callbacks, supply the same callback behavior and inputs.

## Benchmark

Run the collision scenario twice in one process and compare a state hash at
each complete tick, including warmup ticks.

```powershell
.\build\bin\Release\sqphysics2d_benchmark.exe circle_box_bounce --bodies 12 --warmup 2 --ticks 20 --check-determinism
```

Run the executable without a scenario for available scenarios and options.
The reported tick timings measure `World::step()`.

## CMake integration

Install the static library, headers, and CMake package from the Release build.

```powershell
cmake --install build --config Release --prefix build/install
```

Configure your own CMake project with `CMAKE_PREFIX_PATH` pointing to the
install prefix, and link the exported target.

```cmake
find_package(sqphysics2d CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE sqp::sqphysics2d)
```

Include `<sqphysics2d/world.h>` to create a world, submit commands, and capture
or restore snapshots.

## License

- [MIT License](LICENSE)
- [Third-party notices](docs/THIRD_PARTY_NOTICES.md)
