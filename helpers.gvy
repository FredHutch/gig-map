// Function to check for required parameters
def require_param(param_val, param_key) {

    // If the parameter was not provided
    if (!param_val){
        java.util.logging.Logger.getLogger("").info"""
        -----------------
        MISSING PARAMETER: ${param_key}
        -----------------

        Run with --help for more information
        """

        System.exit 1
    }
}

// Function which prints help message text
def help_message(msg, show) {

    // If the user uses the --help flag
    if (show){

        // Show the help message
        java.util.logging.Logger.getLogger("").info"""
        ${msg}
        """.stripIndent()

        // Exit out and do not run anything else
        System.exit 0
    }

}