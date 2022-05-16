from monty.json import MSONable


class _BaseParameters(MSONable):
    def validate(self, structure, sites):
        raise NotImplementedError

    def write(self):
        raise NotImplementedError
