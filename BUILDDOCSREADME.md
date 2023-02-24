## how to build documenation using sphinx
1. cd into docs folder 
2. Open up conf.py and change version number
2. run `sphinx-build -b html . builddir`
    - if you get error that sphinxcontrib-bibtex is missing run pip3 install sphinxcontrib-bibtex
    - if you get error that theme furo is missing then run pip3 install furo
3. The command should complete successfully and produce a directory called builddir, inside of which are the build html files
4. Make a copy of that directory and store temporarily outside of git repo
5. Commit the files
6. git checkout gh-pages 
7. Delete all of the files in repo except .nojekyll and .git
8. Copy all files inside of temporary directory we recently made into the repo. 
9. Commit files
10. git push origin gh-pages
11. Github should automatically update website
12. Visit url and reload page. It should display updated documentation.