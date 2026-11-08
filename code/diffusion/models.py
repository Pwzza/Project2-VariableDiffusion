from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Model:
    name: str = 'smooth'
    contrast: float = 10.0

    def __post_init__(self):
        if self.name not in ('smooth', 'layered') or self.contrast <= 0:
            raise ValueError('Expected smooth/layered and positive contrast')

    def a(self, x, y):
        if self.name == 'smooth':
            return 1.0 + np.asarray(x) + np.zeros_like(y)
        return np.where(np.asarray(x) < .5, 1., self.contrast) + np.zeros_like(y)

    def exact(self, x, y):
        if self.name == 'smooth':
            return np.sin(np.pi*x)*np.sin(np.pi*y)
        c = 2*self.contrast/(1+self.contrast)
        return np.where(np.asarray(x) <= .5, c*np.asarray(x),
                        c*.5+c/self.contrast*(np.asarray(x)-.5)) + np.zeros_like(y)

    def gradient(self, x, y):
        x, y = np.broadcast_arrays(x, y)
        if self.name == 'smooth':
            return np.array([np.pi*np.cos(np.pi*x)*np.sin(np.pi*y),
                             np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)])
        c = 2*self.contrast/(1+self.contrast)
        return np.array([np.where(x < .5, c, c/self.contrast), np.zeros_like(x)])

    def forcing(self, x, y):
        if self.name == 'smooth':
            return (2*np.pi**2*(1+x)*np.sin(np.pi*x)*np.sin(np.pi*y)
                    - np.pi*np.cos(np.pi*x)*np.sin(np.pi*y))
        return np.zeros(np.broadcast_shapes(np.shape(x), np.shape(y)))

    def flux(self, x, y):
        return -self.a(x, y)*self.gradient(x, y)
