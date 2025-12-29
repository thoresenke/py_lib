class MatlabColors:
    def __init__(self):
        # Standard MATLAB colors (R2014b+)
        self._colors = [
            "#0072BD", # 1: Blue
            "#D95319", # 2: Orange
            "#EDB120", # 3: Yellow
            "#7E2F8E", # 4: Purple
            "#77AC30", # 5: Green
            "#4DBEEE", # 6: Light Blue
            "#A2142F"  # 7: Red
        ]

    def __call__(self, idx):
        """
        Get color by 1-based index (MATLAB style).
        Cycles if index > 7.
        """
        if isinstance(idx, int):
            # Adjust for 1-based index and cycle
            return self._colors[(idx - 1) % len(self._colors)]
        raise ValueError("Index must be an integer")

    def c(self, idx):
        """Alias for calling the object directly."""
        return self(idx)

# Create a singleton instance
mCol = MatlabColors()

