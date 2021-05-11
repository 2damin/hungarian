# Hungarian Algorithm
Hungarian Algorithm 

## Build

### Windows

```bash
mkdir build
cd build
cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release ..;cmake --build . --config "Release" -j;
cd ..
```

### Run
```bash
./bin/${CMAKE_BUILD_TYPE}/hungarian.exe
```

## Description

- **MODE** : 0(Mimize total cost), 1(Maximize total cost)
- **N** : the number of workers
- **M** : the number of tasks