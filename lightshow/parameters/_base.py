from abc import ABC, abstractmethod

from monty.json import MSONable


class _BaseParameters(MSONable, ABC):
    @abstractmethod
    def validate(self, structure, sites):
        ...

    @abstractmethod
    def write(self):
        ...

    @property
    @abstractmethod
    def calculation_name(self):
        ...
