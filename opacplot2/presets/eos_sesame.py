#!/usr/bin/python
# -*- coding: utf-8 -*-

SESAME = {3720:
             {'merge': {'filter_dens': lambda x: (x!=0.81)*(x>0),
                        'filter_temps': lambda x: (x>0)},
              'fracs': (1.0,)
              }

         }


