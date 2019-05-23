from itertools import product
for ((x,y),z) in product(
                product(
                  range(2),
                  range(3)),
                range(4)):
    print(x,y,z)

