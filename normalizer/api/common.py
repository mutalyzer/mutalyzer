from flask_restx import abort


def errors(endpoint):
    def dec(*args, **kwargs):
        output = endpoint(*args, **kwargs)
        if output.get("errors"):
            abort(422, "Errors encountered. Check the 'custom' field.", custom=output)
        return output

    return dec
