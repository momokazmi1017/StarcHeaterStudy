"""First-order model of an adaptive-emissivity (STARC-like) radiator in LEO.

FIRST-ORDER PUBLIC-DATA ESTIMATE. NOT VALIDATED AGAINST STARC TEST DATA.
"""

from .params import default_params
from .model import emissivity, stack_coatings, config, simulate, size_radiator, compare
from .economics import radiated_flux, datacenter_area, economics

__all__ = ["default_params", "emissivity", "stack_coatings", "config", "simulate",
           "size_radiator", "compare", "radiated_flux", "datacenter_area", "economics"]
