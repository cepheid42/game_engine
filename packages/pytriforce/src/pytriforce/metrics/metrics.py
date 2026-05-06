from dataclasses import dataclass
from enum import StrEnum

from ..params.header_utils import *


class MetricType(StrEnum):
    ParticleDump = 'MetricType::ParticleDump'
    ParticleDiagnostics = 'MetricType::ParticleDiag'
    ParticleEnergy = 'MetricType::ParticleEnergy'
    FieldDump = 'MetricType::FieldDump'
    FieldEnergy = 'MetricType::FieldEnergy'


class FieldSlice:
    component : str
    starts : tuple
    stops : tuple
    steps : tuple


@dataclass
class Metrics:
    data_path: str = ''
    metrics: tuple = ()
    field_slices: tuple = ()

    def __repr__(self):
        '''

        Check the MDSpan version for best way to define compile time mdspans

        :return:
        '''

        return str(
            section_label('Metrics Parameters') +
            enum_declaration('MetricType', MetricType) +
            '\n' +
            constexpr_declaration('metric_data_path', f'{self.data_path}sv') +
            constexpr_declaration('metric_spec', tuple_to_array(self.metrics), typestr=f'std::array<MetricType, {len(self.metrics)}>')

        )
