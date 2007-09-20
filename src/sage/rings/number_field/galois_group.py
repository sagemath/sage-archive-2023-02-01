

class GaloisGroup:
    def __init__(self, group, number_field):
        self.__group = group
        self.__number_field = number_field

    def __repr__(self):
        return "Galois group %s of the number field %s"%(
            self.__group, self.__number_field)

    def group(self):
        return self.__group

    def order(self):
        return self.__group.order()

    def number_field(self):
        return self.__number_field
