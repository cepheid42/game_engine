from dataclasses import dataclass
from enum import StrEnum

from ..header_utils import *


class MetricType(StrEnum):
    ParticleDump = 'MetricType::ParticleDump'
    ParticleDiagnostics = 'MetricType::ParticleDiag'
    ParticleEnergy = 'MetricType::ParticleEnergy'
    FieldDump = 'MetricType::FieldDump'
    FieldEnergy = 'MetricType::FieldEnergy'


# @dataclass
# class FieldSlice:
#     component : str
#     starts : tuple
#     stops : tuple
#     steps : tuple
#
#     def __repr__(self):
#         extents = (
#             self.stops[0] - self.starts[0],
#             self.stops[1] - self.starts[1],
#             self.stops[2] - self.starts[2],
#         )
#         return (
#             '   FieldSlice{\n'
#             f'      .component = "{self.component}",\n'
#             f'      .starts = {tuple_to_array(self.starts)},\n'
#             f'      .strides = {tuple_to_array(self.steps)},\n'
#             f'      .extents = {tuple_to_array(extents)},\n'
#             '   }'
#         )

@dataclass
class Metrics:
    data_path: str = ''
    metrics: tuple = ()
    # field_slices: tuple = ()

    def __repr__(self):
        return str(
            section_label('Metrics Parameters') +
            enum_declaration('MetricType', MetricType) +
            '\n' +
            constexpr_declaration('metric_data_path', f'{self.data_path}sv') +
            constexpr_declaration('metric_spec', tuple_to_array(self.metrics), typestr=f'std::array<MetricType, {len(self.metrics)}>')

        )
