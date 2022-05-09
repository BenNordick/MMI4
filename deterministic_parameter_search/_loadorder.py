import tellurium
import rrplugins

# Importing this module first ensures that tellurium and rrplugins are loaded in the correct order,
# preventing a double-free crash at exit, even if one file does not need both modules.