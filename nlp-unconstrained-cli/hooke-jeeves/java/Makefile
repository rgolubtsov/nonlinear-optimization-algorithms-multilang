#
# nlp-unconstrained-cli/hooke-jeeves/java/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2025
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# (See the LICENSE file at the top of the source tree.)
#

# Vars from pom.xml (artifactId, version, etc.)
POM_ARTIFACT_ID = hooke-jeeves
POM_VERSION     = 0.1
POM_PACKAGING   = jar

POM_XML             = pom.xml
PATH_TO_APP_JAR     = target
PATH_TO_APP_CLASSES = $(PATH_TO_APP_JAR)/classes
JAVA_PACKAGE_PATH   = optimization/nonlinear/unconstrained/cli
JAR                 = $(POM_ARTIFACT_ID)-$(POM_VERSION).$(POM_PACKAGING)
DEPS                = $(PATH_TO_APP_CLASSES)/$(JAVA_PACKAGE_PATH)/Rosenbrock.class \
                      $(PATH_TO_APP_CLASSES)/$(JAVA_PACKAGE_PATH)/Woods.class \
                      $(PATH_TO_APP_CLASSES)/$(JAVA_PACKAGE_PATH)/Hooke.class \
                      $(PATH_TO_APP_JAR)/$(JAR)
SCRIPT1             = hooke
SCRIPT2S            = -woods
SCRIPT2             = $(SCRIPT1)$(SCRIPT2S)
BIN_DIR             = bin
LIB_DIR             = lib

# Specify flags and other vars here.
# Apache Maven:
MVN         = mvn
MVN_PACKAGE = package
MVN_CLEAN   = clean
MVNFLAGS    = --strict-checksums

# Java:
JAVA       = java
JAVAFLAGS1 = -jar
JAVAFLAGS2 = -cp

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
define SCRIPT1_FULL
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-cli/$(POM_ARTIFACT_ID)/$(JAVA)/$(BIN_DIR)/$(SCRIPT1)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1\n\
# =============================================================================\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r JAVA=$(JAVA)\n\
declare -r JAVAFLAGS=$(JAVAFLAGS1)\n\
declare -r LIB_DIR=$(LIB_DIR)\n\
declare -r JAR=$(JAR)\n\
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
# vim:set nu et ts=4 sw=4:\n\
# =============================================================================
endef

# Woods:
define SCRIPT2_FULL
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-cli/$(POM_ARTIFACT_ID)/$(JAVA)/$(BIN_DIR)/$(SCRIPT2)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1\n\
# =============================================================================\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r JAVA=$(JAVA)\n\
declare -r JAVAFLAGS=$(JAVAFLAGS2)\n\
declare -r LIB_DIR=$(LIB_DIR)\n\
declare -r JAR=$(JAR)\n\
declare -r WOODS_CLASS=optimization.nonlinear.unconstrained.cli.Woods\n\
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
# vim:set nu et ts=4 sw=4:\n\
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
				$(ECHO) $(ECHOFLAGS) "$(SCRIPT1_FULL)" > $(BIN_DIR)/$(SCRIPT1); \
				$(ECHO) $(ECHOFLAGS) "$(SCRIPT2_FULL)" > $(BIN_DIR)/$(SCRIPT2); \
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

# vim:set nu ts=4 sw=4:
