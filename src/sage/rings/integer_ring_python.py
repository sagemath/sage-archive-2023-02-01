def iterator(self):
    yield self(0)
    n = self(1)
    while True:
        yield n
        yield -n
        n += 1

