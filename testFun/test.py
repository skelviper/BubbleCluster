from multiprocessing import Pool

def f(x):
    return x*x

if __name__ == '__main__':
    with Pool(8) as p:
        for i in range(10000000):
            print(p.map(f, [1, 2, 3,4,5,6]))