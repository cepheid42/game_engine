from dataclasses import dataclass
from enum import StrEnum

class MetricType(StrEnum):
    ParticleDump = 'MetricType::ParticleDump'
    ParticleDiagnostics = 'MetricType::ParticleDiag'
    ParticleEnergy = 'MetricType::ParticleEnergy'
    FieldDump = 'MetricType::FieldDump'
    FieldEnergy = 'MetricType::FieldEnergy'


@dataclass
class Metrics:
    data_path: str = ''
    metrics: tuple = ()

    def __repr__(self):
        return ',\n\t'.join(m for m in self.metrics)
