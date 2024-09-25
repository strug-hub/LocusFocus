class PipelineStage():
    """
    An abstract pipeline stage. 
    
    Other stages must inherit from this one and implement the `invoke` and `name` methods.
    """

    def name(self) -> str:
        """
        The unique name of the stage.

        Use kebab-case for the name. eg. "my-example-stage"
        """
        raise NotImplementedError()

    def invoke(self, payload: object) -> object:
        """
        Invoke this pipeline stage on the given payload.

        Return a new or modified payload.
        """
        raise NotImplementedError()
