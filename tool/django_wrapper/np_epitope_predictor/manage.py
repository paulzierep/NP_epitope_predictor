#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys

############################
#The modules need to be loaded in the manage.py for ppickle to work
############################

from np_epitope_predictor.settings import PREDICTOR_PATH

import sys
sys.path.insert(0, PREDICTOR_PATH)

from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor

############################
#Noral django setup
############################

def main():
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'np_epitope_predictor.settings')
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    main()
