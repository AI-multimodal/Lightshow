from abc import ABC, abstractmethod


class _BaseParameters(ABC):
    @abstractmethod
    def write(self, target_directory, kwargs):
        ...
