from abc import ABC, abstractmethod

from monty.json import MSONable


class _BaseParameters(MSONable, ABC):
    @abstractmethod
    def write(self, target_directory, kwargs):
        ...

    @property
    @abstractmethod
    def calculation_name(self):
        ...
