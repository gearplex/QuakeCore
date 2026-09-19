"""Exhaustive numerical completion labels; no label implies physical collapse."""


def completion_status(quakecore_status, opensees_status, quakecore_steps, opensees_steps):
    allowed = {"completed", "numerical_noncompletion"}
    if quakecore_status not in allowed or opensees_status not in allowed:
        raise ValueError("unknown solver completion status")
    for steps in (quakecore_steps, opensees_steps):
        if not isinstance(steps, int) or isinstance(steps, bool) or steps < 0:
            raise ValueError("committed nominal step counts must be nonnegative integers")
    qdone = quakecore_status == "completed"
    odone = opensees_status == "completed"
    both_stopped = not qdone and not odone
    result = {
        "both_completed": qdone and odone,
        "shared_numerical_noncompletion_same_step": both_stopped and quakecore_steps == opensees_steps,
        "dual_numerical_noncompletion_different_steps": both_stopped and quakecore_steps != opensees_steps,
        "opensees_reference_limited": qdone and not odone,
        "quakecore_limited": not qdone and odone,
    }
    assert sum(result.values()) == 1
    return result
