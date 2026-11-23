# SConstruct

import glob
import os

# Options

vars = Variables(None, ARGUMENTS)

vars.Add(EnumVariable(
    "target",
    "Build target",
    "debug",
    allowed_values=("debug", "release"),
))

vars.Add(EnumVariable(
    "precision",
    "Default real_t precision",
    "float",
    allowed_values=("float", "double"),
))

vars.Add(BoolVariable(
    "tests",
    "Build tests",
    True,
))

vars.Add(BoolVariable(
    "examples",
    "Build examples",
    True,
))

vars.Add(BoolVariable(
    "validation",
    "Runtime debug validations",
    True,
))

vars.Add(BoolVariable(
    "shared",
    "Build shared library",
    False,
))

vars.Add(PathVariable(
    "prefix",
    "Install prefix",
    "/usr/local",
    PathVariable.PathAccept,
))

env = Environment(
    variables=vars,
    COMPILATIONDB_USE_ABSPATH=True,
)

env.Tool("compilation_db")

Help(vars.GenerateHelpText(env))

# Basic config

target = env["target"]
precision = env["precision"]
validation = env["validation"]

project_name = "sqphysics2d"

build_dir = os.path.join("build", target)
lib_dir = os.path.join(build_dir, "lib")
bin_dir = os.path.join(build_dir, "bin")

env.Append(CPPPATH=[
    ".",
    "include",
])

env.Append(CPPDEFINES=[
    "SQP_BUILD",
])

if precision == "double":
    env.Append(CPPDEFINES=["SQP_REAL_DOUBLE"])
else:
    env.Append(CPPDEFINES=["SQP_REAL_FLOAT"])

if validation:
    env.Append(CPPDEFINES=["SQP_DEBUG_VALIDATION"])

# Compiler flags

is_msvc = env.get("MSVC_VERSION") is not None

if is_msvc:
    env.Append(CXXFLAGS=[
        "/std:c++20",
        "/W4",
        "/EHsc",
    ])

    if target == "debug":
        env.Append(CXXFLAGS=[
            "/Zi",
            "/Od",
            "/DDEBUG",
        ])
        env.Append(LINKFLAGS=["/DEBUG"])
    else:
        env.Append(CXXFLAGS=[
            "/O2",
            "/DNDEBUG",
        ])
else:
    env.Append(CXXFLAGS=[
        "-std=c++20",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wconversion",
        "-Wshadow",
    ])

    if target == "debug":
        env.Append(CXXFLAGS=[
            "-g",
            "-O0",
            "-DDEBUG",
        ])
    else:
        env.Append(CXXFLAGS=[
            "-O3",
            "-DNDEBUG",
        ])

# Source collection

def glob_sources_recursive(path):
    return sorted(glob.glob(os.path.join(path, "**", "*.cpp"), recursive=True))

sources = glob_sources_recursive("src")

# Library

lib_env = env.Clone()
lib_env.Append(CPPPATH=["src"])

if env["shared"]:
    sqp_lib = lib_env.SharedLibrary(
        target=os.path.join(lib_dir, project_name),
        source=sources,
    )
else:
    sqp_lib = lib_env.StaticLibrary(
        target=os.path.join(lib_dir, project_name),
        source=sources,
    )

Default(sqp_lib)

# Examples
example_program = None
if env["examples"]:
    example_env = env.Clone()
    example_env.Append(LIBPATH=[lib_dir])
    example_env.Append(LIBS=[project_name])
    example_program = example_env.Program(
        target=os.path.join(bin_dir, "sqphysics2d_benchmark"),
        source=glob_sources_recursive(os.path.join("examples", "benchmark")),
    )
    Depends(example_program, sqp_lib)
    Default(example_program)

# Tests

test_program = None

if env["tests"] and os.path.isdir("tests"):
    test_env = env.Clone()
    test_env.Append(CPPPATH=["tests", "src", "extern/catch2/extras"])
    test_env.Append(LIBPATH=[lib_dir])
    test_env.Append(LIBS=[project_name])

    test_sources = glob_sources_recursive("tests")
    test_sources.append("extern/catch2/extras/catch_amalgamated.cpp")

    if test_sources:
        test_program = test_env.Program(
            target=os.path.join(bin_dir, "sqphysics2d_tests"),
            source=test_sources,
        )

        Depends(test_program, sqp_lib)
        Default(test_program)

# Install

install_include_dir = os.path.join(env["prefix"], "include")
public_include_dir = os.path.join(install_include_dir, "sqphysics2d")
install_headers = [
    env.InstallAs(os.path.join(public_include_dir, name), os.path.join("include/sqphysics2d", name))
    for name in (
        "build_config.h",
        "settings.h",
        "material.h",
        "shape.h",
        "body.h",
        "joint.h",
        "world.h",
        "events.h",
    )
]
detail_include_dir = os.path.join(public_include_dir, "detail")
install_headers += [
    env.InstallAs(os.path.join(detail_include_dir, name), os.path.join("include/sqphysics2d/detail", name))
    for name in (
        "types.h",
        "math.h",
        "contact.h",
        "contact_solver.h",
        "constraint_graph.h",
        "data_storage.h",
        "core_world.h",
        "broad_phase.h",
        "shape_query_view.h",
    )
]

install_lib = env.Install(
    os.path.join(env["prefix"], "lib"),
    sqp_lib,
)

env.Alias("install", install_headers + install_lib)

# Convenience aliases

env.Alias("lib", sqp_lib)

if test_program is not None:
    env.Alias("test", test_program)
if example_program is not None:
    env.Alias("examples", example_program)

env.Alias("clean-build", [], Delete("build"))

compdb = env.CompilationDatabase("compile_commands.json")
Default(compdb)
env.Alias("compdb", compdb)
