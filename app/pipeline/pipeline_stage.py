from .payload import Payload

class PipelineStage():
    """
    An abstract pipeline stage. Other stages should inherit
    from this one, and implement the `invoke` method.
    """

    def invoke(self, payload: Payload) -> Payload:
        """
        Invoke this pipeline stage on the given payload.

        Return a new or modified payload.
        """
        raise NotImplementedError()
