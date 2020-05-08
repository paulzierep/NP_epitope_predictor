from django.urls import path

from . import views

urlpatterns = [
    path(r'', views.index, name='index'),
    path(r'index', views.index, name='index'),
    path(r'help', views.help, name='help'),
    path(r'results/<int:file_id>', views.results, name='results'),
]