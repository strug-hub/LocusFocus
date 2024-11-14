class LocusFocusError(Exception):
    """
    Generic error class for LocusFocus.
    """

    status_code = 500

    def __init__(self, message, status_code=None, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv["message"] = self.message
        return rv

    def __str__(self):
        return f"{self.status_code} {self.message}"


class InvalidUsage(LocusFocusError):
    """
    Exceptions raised due to user errors (eg. invalid input).
    """

    status_code = 400


class ServerError(LocusFocusError):
    """
    Exceptions raised due to server errors 
    (eg. database connection issues, missing files, etc.).
    """

    status_code = 500
