
from django.shortcuts import render, redirect
from django.http import HttpResponse

import pickle

from .forms import SMILES_input

from np_epitope_predictor.settings import PREDICTOR_PATH

import sys
sys.path.insert(0, PREDICTOR_PATH)

from run_prediction_wrapper import prediction2django

import os
#from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
# sys.modules['Molecule_Group_Classifier'] = Molecule_Group_Classifier

from np_epitope_predictor.settings import MEDIA_ROOT

def index(request):

    if not request.method == 'POST':

        form = SMILES_input()
        return render(request, 'index.html', {'form':form})
    else:
        form = SMILES_input(request.POST) #create a form object
        if form.is_valid():
            smiles = form.cleaned_data['smiles']
            results = prediction2django(smiles)

            output_path = os.path.join(MEDIA_ROOT, 'results.pickle')

            with open(output_path, 'wb') as handle:
                pickle.dump(results, handle)

            return redirect('results')
            # return render(request, 'results.html', {'results':results})


def results(request):


    output_path = os.path.join(MEDIA_ROOT, 'results.pickle')

    with open(output_path, 'rb') as handle:
        results = pickle.load(handle)

    return render(request, 'results.html', {'results':results})