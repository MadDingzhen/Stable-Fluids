# Stable-Fluids
This is the implementation of *Stable Fluids* by Jos Stam at SIGGRAPH 1999.

FFT method with periodic boundary condition is one way to solve the equation of diffuse and projection faster.


## Build
**Windows is not supported.**

Libraries you need are
- FFTW (version 3)
- GLM
- GLEW
- GLFW

You can install the required libraries by running the following command.

```shell
$ sudo apt-get update && sudo apt-get install libxi-dev libxinerama-dev libxcursor-dev libglm-dev libglfw3-dev libfftw3-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev libglew-dev
```


You can build and execute following command.

```shell
$ git clone https://github.com/daichi-ishida/Stable-Fluids.git
$ make
$ ./bin/main
```

Density is continuously added to the initial area by default.

If you want to set density only once at the beginning, 
add `-o` option that represents for "once".

```shell
$ ./bin/main -o
```

## Screenshot

![Screenshot](output/Screenshot.png)

Each of red lines shows the direction of the grid's velocity and 
brightness represents density at the position.