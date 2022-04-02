# Getting started

This guide gives a brief overview of how this library is organized, and how you can create VPC scenarios to simulate.

## Structure

The library is organized as follows:

```sh
├── core/           # VPC essential code files library
│   ├── anim/       # Animation functions go here
│   ├── gp/         # Functions related to GP go here
│   ├── rbm/        # Rigid Body Motion functions
│   ├── simulink/   # Simulink blocks for VPC are defined here
│   ├── trafo/      # Functions related to transformations
│   └── vision/     # Functions related to vision/cameras
├── gpml/           # Gaussian Process library (GPML)
├── templates/      # VMO, VPC, GP... examplary simulink template files 
└── tests/          # Unit tests to check for failures during library development
```

## How to use the library

Before using any files in the library, be sure to run the `startup.m` script in the root directory of the repository.
This script will initialize and load the VPC block library into the Simulink Block library window.
It will also load all scripts into your `PATH` so that the functions within the library can find each other.

Since the library uses internally the [GPML library](http://www.gaussianprocess.org/gpml/code/matlab/doc/) for Gaussian Process Regression, be sure that the submodule has been cloned to your local machine as well. You can update all submodules by the following command:

```bash
git submodule update --recursive --remote
```

> If it is your first time cloning the repository, use this command instead:
>
> ```bash
> git submodule update --init --recursive
> ```

### Example files

Inside the `templates/` folder you will find Simulink templates for a simple Visual Motion Observer (`VMO.slx`) and how to use it in a Visual Pursuit Control scenario (`VPC.slx`). Simply run them and be sure to understand them so that you can use them for your own purpose.

### Library blocks

The VPC library blocks are defined within:

```sh
core/simulink/
├── simlib_vpc.slx          # ROOT Library header that links all sub-libraries below
├── simlib_vpc_basic.slx    # Basic VPC-blocks (VMO, Virtual Camera, RBM Model, ...)
├── simlib_vpc_coremath.slx # VPC core blocks related to transformations, vector operations, ...
└── simlib_vpc_gp.slx       # GP block library
```

If you need to extend the VPC library, you can do it either there, or define your own library and link to the VPC library so that they show up as individual blocks in your library (recommended approach).

## FAQ

> The VPC library does not show up in the Simulink block library window!

Ensure that you have run the `startup.m` script inside the repository's root folder. Then, refresh the Simulink block library window by pressing `F5` or `right-click -> Refresh Window`.

> I opened the template Simlink files, but some blocks have a "`?`" written on them, and I get a lot of warnings.

Run the `startup.m` script again and re-open the template file.

> I tried all of the solutions above in my own code directory but nothing worked.

Check your paths. If you type `which rbm` in the MATLAB console and nothing appears, the VPC library has not been fully loaded. Try going to the VPC library root folder within MATLAB and then run the `startup.m` script within there (rather than within your project). In any way, it must be a `PATH` error.
