from unittest.mock import MagicMock

import pytest
from pymongo.errors import AutoReconnect

from app.utils.gtex_db.real import (
    RealGTExDatabase,
    _RECONNECT_ATTEMPTS,
    _retry_on_disconnect,
)


def test_retry_on_disconnect_recovers_before_attempts_exhausted(monkeypatch):
    """A transient AutoReconnect is retried and the call eventually succeeds."""
    monkeypatch.setattr("app.utils.gtex_db.real.time.sleep", lambda _: None)

    calls = {"n": 0}

    @_retry_on_disconnect
    def flaky():
        calls["n"] += 1
        if calls["n"] < _RECONNECT_ATTEMPTS:
            raise AutoReconnect("connection lost")
        return "reconnected"

    assert flaky() == "reconnected"
    assert calls["n"] == _RECONNECT_ATTEMPTS


def test_retry_on_disconnect_raises_after_exhausting_attempts(monkeypatch):
    """If MongoDB never comes back, the original error still propagates."""
    monkeypatch.setattr("app.utils.gtex_db.real.time.sleep", lambda _: None)

    @_retry_on_disconnect
    def always_down():
        raise AutoReconnect("still down")

    with pytest.raises(AutoReconnect):
        always_down()


def test_list_tissues_retries_through_real_gtex_database(monkeypatch):
    """RealGTExDatabase.list_tissues itself retries on AutoReconnect."""
    monkeypatch.setattr("app.utils.gtex_db.real.time.sleep", lambda _: None)

    fake_db = MagicMock()
    fake_db.list_collection_names.side_effect = [
        AutoReconnect("blip"),
        ["Liver", "variant_table"],
    ]
    client = {"GTEx_V8": fake_db}

    db = RealGTExDatabase(client)
    assert db.list_tissues("V8") == ["Liver"]
    assert fake_db.list_collection_names.call_count == 2
