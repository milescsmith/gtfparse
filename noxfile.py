import nox


@nox.session
def run_tests(session):
    session.install(".", "pytest")
    session.run("pytest", "-v", "tests")
