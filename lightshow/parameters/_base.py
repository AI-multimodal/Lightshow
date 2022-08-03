from abc import ABC, abstractmethod


class _BaseParameters(ABC):
    @abstractmethod
    def write(self, target_directory, kwargs):
        ...

    @property
    def name(self):
        return self._name
