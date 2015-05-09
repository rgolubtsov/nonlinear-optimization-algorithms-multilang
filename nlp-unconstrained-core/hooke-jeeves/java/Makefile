# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/java/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Copyright (C) 2015 Radislav (Radic) Golubtsov

# Vars from pom.xml (artifactId, version, etc.)
POM_ARTIFACT_ID = hooke-jeeves
POM_VERSION     = 0.1
POM_PACKAGING   = jar

POM_XML             = pom.xml
PATH_TO_APP_JAR     = target
PATH_TO_APP_CLASSES = $(PATH_TO_APP_JAR)/classes
JAVA_PACKAGE_PATH   = optimization/nonlinear/unconstrained/core
JAR                 = $(POM_ARTIFACT_ID)-$(POM_VERSION).$(POM_PACKAGING)
DEPS                = $(PATH_TO_APP_CLASSES)/$(JAVA_PACKAGE_PATH)/Rosenbrock.class \
                      $(PATH_TO_APP_CLASSES)/$(JAVA_PACKAGE_PATH)/Woods.class \
                      $(PATH_TO_APP_JAR)/$(JAR)

# Specify flags and other vars here.
# Apache Maven:
MVN         = mvn
MVN_PACKAGE = package
MVN_CLEAN   = clean
MVNFLAGS    = --strict-checksums

ECHO         = echo
EXIT         = exit
EXIT_FAILURE = 1
EXIT_SUCCESS = 0

# Making the target.
$(DEPS):
	if [ ! -e "$(POM_XML)" ]; then \
		$(ECHO) "Error: The main build scenario ($(POM_XML)) does not exist. Exiting..."; \
		$(EXIT) $(EXIT_FAILURE); \
	fi

	$(MVN) $(MVNFLAGS) $(MVN_PACKAGE)

.PHONY: all clean

all: $(DEPS)

clean:
	$(MVN) $(MVN_CLEAN)

	$(RM) $(DEPS)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
