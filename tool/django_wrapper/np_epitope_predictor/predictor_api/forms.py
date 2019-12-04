from django import forms

class SMILES_input(forms.Form):
	smiles = forms.CharField(label='SMILES', 
	initial = 'C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@@H](CON)NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@H](O)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(O)=O',
	required = True,
	help_text='Insert SMILES of target molecule.',
	widget=forms.Textarea(attrs={'cols': 100, 'rows': 4})
	)