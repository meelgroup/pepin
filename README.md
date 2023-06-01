# Pepin

Pepin is a DNF streaming counting tool that gives the number of approximate points in a volume given a set of n-dimensional boxes that may intersect. Think of it as a tool that takes a set of cubes in space and approximates the total volume of all cubes. The related research paper is "Engineering an Efficient Approximate DNF-Counter".

It takes in a set of cubes from a DNF file, such as this, that has 10 dimensions, and 2 cubes:

```
$ cat myfile.dnf
p dnf 20 2
1 2 3 0
1 -5 0
```

And outputs a probabilistically approximate count.

You can build&run it with:

```
mkdir build
cd build
cmake ..
make
./pepin --epsilon 0.15 --delta 0.1 myfile.dnf
[...]
Approx num points: 339456
```

Notice that the cube `1 2 3` has 2**17=131072 points, and `1 -5` has 2**18=262144, so a total of 393216. But they overlap, `1 2 3 -5` is counted twice. So the exact number is: 327680. Hence, we over-approximated a bit here. The error is 1.0-339456/327680=.0359, so about 3.6%. This is well below the advertised 15% error allowed (i.e. epsilon 0.15).


## Current Limitations

The following will likely be lifted in 1-2 weeks:
* `make install` is not yet working, will fix it next week.
* Number of variables must be divisible by 4
* [May take more time] Can only return estimate after the exact number of clauses have been passed in as promised

The last limitation means you MUST have the `p dnf VARS CLS` header correct in your DNF file or the tool will NOT work.

## Library Use


The header file, `pepin.h` is made to be used as a library. Use as:
```
#include <pepin.h>
#include <iostream>
#include <iomanip>
using namespace PepinNS;

int main() {
  Pepin pepin;
  pepin.new_vars(20);
  pepin.set_n_cls(2);

  typedef vector<Lit> CL;
  pepin.add_clause(CL{itol(1), itol(2), itol(3)});
  pepin.add_clause(CL{itol(1), itol(-5)});
  auto weigh_num_sols = dnfs->get_low_prec_appx_weighted_sol();
  std::cout << "Solution: " << std::scientific << std::setprecision(10)
    << *weigh_num_sols << std::endl;

  return 0;
}
```

The library is clean -- it cleans up after itself, you can have more than one in memory, you can even run them in parallel, if you wish.

## Fuzzing

You can fuzz by building DNFKLM from [here](https://gitlab.com/Shrotri/DNF_Counting/), putting the resulting `DNFKLM` binary into `build/`, and then:

```
cd build
ln -s ../scripts/* .
./build_norm.sh
mkdir -f tests
./fuzz_test.py
```

You can check `fuzz_test.py` to adjust etc.

## License
MIT license all the way through
