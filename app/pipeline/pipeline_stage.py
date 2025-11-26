class PipelineStage:
    """
    An abstract pipeline stage.

    Other stages must inherit from this one.
    """

    def name(self) -> str:
        """
        The unique name of the stage.

        Use kebab-case for the name. eg. "my-example-stage"
        """
        raise NotImplementedError()

    def description(self) -> str:
        """
        A human-readable description of the stage.

        By default, the same as name but as a sentence.
        """
        return self.name().replace("-", " ").title()

    def invoke(self, payload: object) -> object:
        """
        Invoke this pipeline stage on the given payload.

        Return a new or modified payload.
        """
        raise NotImplementedError()
