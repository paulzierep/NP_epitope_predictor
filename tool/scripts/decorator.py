import pandas as pd

def as_html(param):
	def decorator(func):
		def html_wrapper(*args, **kwargs):
			df = func(*args, **kwargs)
			df = df.to_html()
			return(df)

		return(html_wrapper)
	return(decorator)

@as_html(2)
def create_df(val):

	df = pd.DataFrame()
	df.loc[0,0] = val

	return(df)

df = create_df(1)
print(df)
