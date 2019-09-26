from ..layout import ValRepr


def test_init_ValRepr():
    valrep = ValRepr('erg')
    assert valrep.valrepr['cutoff_MK'] == '$E_\mathrm{Kp,flare,}$$_\mathrm{min}\;[$erg$]$'
    assert valrep.valrepr['truncated_2'] == '$tr_2$'
    assert isinstance(valrep.title, str)
    valrep = ValRepr('s')
    assert valrep.valrepr['cutoff_MK'] == '$ED$$_\\mathrm{min}\\;[$s$]$'
    assert valrep.valrepr['truncated_2'] == '$tr_2$'
    assert isinstance(valrep.title, str)