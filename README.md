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

Notice that the cube `1 2 3` has 2**17 points, so 131072 points, and `1 -5` has 2**18, so 262144, so a total of 393216. But they overlap, `1 2 3 -5` is counted twice. So the exact number is: 327680. Hence, we over-approximated a bit.

## Library Use

The header file, `pepin.h` is made to be used as a library. Use as:
```
#include <pepin.h>
using namespace PepinNS;

int main() {
  Pepin pepin;
  pepin.new_vars(20);
  pepin.set_n_cls(2);

  vector<Lit> cl = {itol(1), itol(2), itol(3)};
  pepin.add_clause(cl);
    
}
```


## License
MIT license all the way through
