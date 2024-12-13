

class ProcessError(Exception):
    pass


def chain_exception_message(message, exception):

    if exception is not None:
        exception_message = str(exception)
        if exception_message and not exception_message.isspace():
            return (message[:-1] if message.endswith(".") else message) + ": " + \
                   (exception_message if exception_message.endswith(".") else exception_message + ".")
    return message if message.endswith(".") else message + "."


def illegal_location_reached():

    raise AssertionError("Program control flow reached a location that should not be reachable.")
