#/usr/bin/env python

class Config(object):

    PROGRAM_OUT =  '/home/licon/uni-koeln/tr32/stats/size_distributions'
    TEST_OUT = '/home/licon/uni-koeln/tr32/homology/homology/components/test/data'

class LES(Config):

    PROGRAM_OUT = '%s/LES' % Config.PROGRAM_OUT


class DALES(Config):

    PROGRAM_OUT = '%s/DALES' % Config.PROGRAM_OUT


class DNS(Config):

    PROGRAM_OUT = '%s/DNS' % Config.PROGRAM_OUT