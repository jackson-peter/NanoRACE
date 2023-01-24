# ADD ALL THE FLEPSEQ SCRIPT (FROM THE GITHUB REPO) TO THE PATH
FLEPSEQ_PATH = str(Path(workflow.basedir) / "scripts")
os.environ["PATH"] += os.pathsep + FLEPSEQ_PATH