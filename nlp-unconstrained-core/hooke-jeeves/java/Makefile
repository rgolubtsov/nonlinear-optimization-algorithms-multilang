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
BIN_DIR             = bin
LIB_DIR             = lib

# Specify flags and other vars here.
# Apache Maven:
MVN         = mvn
MVN_PACKAGE = package
MVN_CLEAN   = clean
MVNFLAGS    = --strict-checksums

ECHO         = echo
ECHOFLAGS    = -e
EXIT         = exit
EXIT_FAILURE = 1
EXIT_SUCCESS = 0
MKDIR        = mkdir
SED          = sed
SEDFLAGS     = -i -e
SEDREGEX     = "s/^\ //g;s/'//g"
CHMOD        = chmod
CHMODFLAGS   = -v
CHMODMODE    = 755
CP           = cp
CPFLAGS      = -v
RMFLAGS      = -vR

# Startup scripts for Rosenbrock and Woods test problems like here documents.
# Rosenbrock:
HOOKE = hooke
define HOOKE_SCRIPT
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-core/hooke-jeeves/java/bin/hooke\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
# Nonlinear programming algorithms as the (un-)constrained minimization\n\
# problems with the focus on their numerical expression using various\n\
# programming languages.\n\
#\n\
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.\n\
# =============================================================================\n\
# Copyright (C) 2015 Radislav (Radic) Golubtsov\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r JAVA=java\n\
declare -r JAVAFLAGS=-jar\n\
declare -r LIB_DIR=lib\n\
declare -r JAR=hooke-jeeves-0.1.jar\n\
\n\
$$'{JAVA}' $$'{JAVAFLAGS}' $$'{LIB_DIR}'/$$'{JAR}'\n\
\n\
if [ $$'?' -eq $$'{EXIT_FAILURE}' ]; then\n\
    exit $$'{EXIT_FAILURE}'\n\
fi\n\
\n\
exit $$'{EXIT_SUCCESS}'\n\
\n\
# =============================================================================\n\
# vim:set nu:et:ts=4:sw=4:\n\
# =============================================================================
endef

# Woods:
HOOKE_WOODS = hooke-woods
define HOOKE_WOODS_SCRIPT
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-core/hooke-jeeves/java/bin/hooke-woods\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
# Nonlinear programming algorithms as the (un-)constrained minimization\n\
# problems with the focus on their numerical expression using various\n\
# programming languages.\n\
#\n\
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.\n\
# =============================================================================\n\
# Copyright (C) 2015 Radislav (Radic) Golubtsov\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r JAVA=java\n\
declare -r JAVAFLAGS=-cp\n\
declare -r LIB_DIR=lib\n\
declare -r JAR=hooke-jeeves-0.1.jar\n\
declare -r WOODS_CLASS=optimization.nonlinear.unconstrained.core.Woods\n\
\n\
$$'{JAVA}' $$'{JAVAFLAGS}' $$'{LIB_DIR}'/$$'{JAR}' $$'{WOODS_CLASS}'\n\
\n\
if [ $$'?' -eq $$'{EXIT_FAILURE}' ]; then\n\
    exit $$'{EXIT_FAILURE}'\n\
fi\n\
\n\
exit $$'{EXIT_SUCCESS}'\n\
\n\
# =============================================================================\n\
# vim:set nu:et:ts=4:sw=4:\n\
# =============================================================================
endef

# Making the target.
$(DEPS) $(BIN_DIR) $(LIB_DIR):
	if [ ! -e "$(POM_XML)" ]; then \
		$(ECHO) "Error: The main build scenario ($(POM_XML)) does not exist. Exiting..."; \
		$(EXIT) $(EXIT_FAILURE); \
	fi

	$(MVN) $(MVNFLAGS) $(MVN_PACKAGE)

	if [ ! -d "$(BIN_DIR)" ]; then \
		$(MKDIR) $(BIN_DIR); \
		if [ ! -d "$(LIB_DIR)" ]; then \
			$(MKDIR) $(LIB_DIR); \
		fi; \
		if [ ! -d "$(PATH_TO_APP_JAR)" ]; then \
			$(ECHO) "Error: There is no such directory ($(PATH_TO_APP_JAR)) - from where the app JAR should be taken. Exiting..."; \
			$(EXIT) $(EXIT_FAILURE); \
		else \
			if [ ! -e "$(PATH_TO_APP_JAR)/$(JAR)" ]; then \
				$(ECHO) "Error: The app JAR ($(JAR)) is not found. Exiting..."; \
				$(EXIT) $(EXIT_FAILURE); \
			else \
				$(ECHO) $(ECHOFLAGS) "$(HOOKE_SCRIPT)" > $(BIN_DIR)/$(HOOKE); \
				$(ECHO) $(ECHOFLAGS) "$(HOOKE_WOODS_SCRIPT)" > $(BIN_DIR)/$(HOOKE_WOODS); \
				$(SED) $(SEDFLAGS) $(SEDREGEX) $(BIN_DIR)/*; \
				$(CHMOD) $(CHMODFLAGS) $(CHMODMODE) $(BIN_DIR)/*; \
				$(CP) $(CPFLAGS) $(PATH_TO_APP_JAR)/$(JAR) $(LIB_DIR); \
			fi \
		fi \
	fi

.PHONY: all clean

all: $(DEPS) $(BIN_DIR) $(LIB_DIR)

clean:
	$(MVN) $(MVN_CLEAN)

	$(RM) $(RMFLAGS) $(DEPS) $(BIN_DIR) $(LIB_DIR)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
