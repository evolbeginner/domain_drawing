#! /bin/env python

from biograpy import Panel, tracks, features
from Bio import SeqFeature

panel = Panel(fig_width=900)

test_track = tracks.BaseTrack(name = 'domains', sort_by = None, cm = 'Set3')

feat = SeqFeature.SeqFeature()
feat.location = SeqFeature.FeatureLocation(100, 300)

'''
#dfeat = features.DomainFeature([feat], name = 'test domain 1', seq_line = (1, 766))
dfeat = features.DomainFeature([feat], name = 'test domain 1')
test_track.append(dfeat)

dfeat = features.DomainFeature([feat],name = 'test domain 1', height = 1.5, seq_line = (1, 766))
test_track.append(dfeat)

dfeat = features.DomainFeature([feat],name = 'test domain 1', height = 1.5, boxstyle = 'round4, rounding_size=1.4',  seq_line = (1, 766))
test_track.append(dfeat)

dfeat = features.DomainFeature([feat],name = 'test domain 1',height = 1.5, boxstyle = 'larrow, pad = 0',  seq_line = (1, 766))
test_track.append(dfeat)

dfeat = features.DomainFeature([feat], name = 'test domain 1',height = 1.5, boxstyle = 'round ',  seq_line = (1, 766))
test_track.append(dfeat)
'''

feat2 = SeqFeature.SeqFeature()
feat2.location = SeqFeature.FeatureLocation(500, 700)

#dfeat = features.DomainFeature([feat,feat2], name = 'test domain 1',height = 1.5, boxstyle = 'roundtooth, pad = 0.1, tooth_size=1.2',  seq_line = (1, 766), alpha = 0.7)
#test_track.append(dfeat)

#dfeat = features.DomainFeature([feat,feat2],name = ['test domain 1', 'test domain 2'], height = 1.5, boxstyle = ['sawtooth, tooth_size=1.4',  'rarrow, pad = 0.1'],  seq_line = (1, 766), ec = 'None')
#test_track.append(dfeat)

dfeat = features.DomainFeature([feat,feat2],name = ['', 'test domain 2'],height = 1.5, boxstyle = ['round4, rounding_size=1.4',  'larrow, pad = 0.1'],   ec = 'None', fc = ['y', 'c'], color_by_cm = False)
test_track.append(dfeat)

panel.add_track(test_track)
panel.save('domain_test.png')

