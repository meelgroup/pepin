# Pepin

Pepin is a DNF streaming counting tool that gives the number of approximate
points in a volume given a set of n-dimensional boxes that may intersect. Think
of it as a tool that takes a set of cubes in space and approximates the total
volume of all cubes. The related research paper, published at IJCAI 2023 is
[Engineering an Efficient Approximate
DNF-Counter](https://www.ijcai.org/proceedings/2023/226).

The tool takes in a set of cubes from a DNF file, such as this, that has 10
dimensions, and 2 cubes:
```plain
$ cat myfile.dnf
p dnf 20 2
1 2 3 0
1 -5 0
```

And outputs a probabilistically approximate count.

## Building
We suggest downloading a [released binary](https://github.com/meelgroup/pepin/releases).
You can also use it via Nix: simply [install Nix](https://nixos.org/download/) and then:
```shell
nix shell github:meelgroup/pepin
```
You will then have the `pepin` binary available in your path.

If you want to build it from source, you can do so with the following commands:
```plain
sudo apt-get install libgmp-dev zlib1g-dev cmake build-essential git
git clone https://github.com/meelgroup/pepin
mkdir build
cd build
cmake ..
make
```

## Usage
Run it via:
```plain
./pepin --epsilon 0.15 --delta 0.1 myfile.dnf
[...]
c [dnfs] Low-precision approx num points: 348672
```

Notice that the cube `1 2 3` has `2**17=131072` points, and `1 -5` has
`2**18=262144`, so a total of 393216. But they overlap, `1 2 3 -5` is counted
twice. So the exact number is: 327680. Hence, we over-approximated a bit here.
The error is `1.0-348672/327680=-0.064`, so about 6.4%. This is well below the
advertised 15% error allowed (i.e. epsilon 0.15).

## Weighted Counting
Pepin is a streaming, weighted approximate model. Because it works with data
streams, you must declare the weights of all literals before you start the
stream. Hence, the weights need to be declared at the top of the DNF file. Here
is an example:
```plain
$ cat myfile-weighted.dnf
p dnf 20 2
w 1 2/3
w 2 3/4
1 2 3 0
1 -5 0
```

The weights of variables 1 and 2 are now declared to be 2/3 and 3/4,
respectively. The rest of the variables will have weight 1/2 by default.

## Current Limitations
Currently, the algorithm can only return estimate after the exact number of
clauses have been passed in as promised. This means you must have the `p dnf
VARS CLS` header correct in your DNF file or the tool will not work.

## Library Use
You can install the library with `sudo make install`. Then, the installed
header file `pepin/pepin.h` can be used as a library:
```cpp
#include <pepin/pepin.h>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace PepinNS;

int main() {
  Pepin pepin;
  pepin.new_vars(20);
  pepin.set_n_cls(2);

  typedef std::vector<Lit> CL;
  pepin.add_clause(CL{itol(1), itol(2), itol(3)});
  pepin.add_clause(CL{itol(1), itol(-5)});
  auto weigh_nsols = pepin.get_low_prec_appx_weighted_sol();
  std::cout << "Solution: " << std::scientific << *weigh_nsols << "std::endl;

  return 0;
}
```

The library is clean -- it cleans up after itself, you can have more than one
in memory, you can even run them in parallel, if you wish.

## Fuzzing
You can fuzz by building DNFKLM from
[here](https://gitlab.com/Shrotri/DNF_Counting/), putting the resulting
`DNFKLM` binary into `build/`, and then:
```bash
cd build
ln -s ../scripts/* .
./build_norm.sh
mkdir -f tests
./fuzz_test.py
```

You can check `fuzz_test.py` to adjust etc.

## License
MIT license all the way through
