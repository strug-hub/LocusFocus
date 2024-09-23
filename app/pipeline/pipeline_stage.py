class PipelineStage():
    """
    An abstract pipeline stage. 
    
    Other stages must inherit from this one and implement the `invoke` method.
    """

    def name(self) -> str:
        """
        The name of the stage.
        """
        raise NotImplementedError()

    def invoke(self, payload: object) -> object:
        """
        Invoke this pipeline stage on the given payload.

        Return a new or modified payload.
        """
        raise NotImplementedError()
