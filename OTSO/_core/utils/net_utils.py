from __future__ import annotations

import os
import socket
from urllib.parse import urlsplit

import requests
import urllib3.util.connection as _u3conn

CONNECT_TIMEOUT = 6        # seconds per address; the read timeout is the caller's `timeout`
IPV6_PROBE_TIMEOUT = 2     # seconds, one-off check per host

_PROBED: set[str] = set()


def prefer_ipv4() -> None:
    """Resolve hosts to IPv4 only (affects requests/urllib3 for the rest of the process)."""
    _u3conn.allowed_gai_family = lambda: socket.AF_INET


def ipv4_only() -> bool:
    return _u3conn.allowed_gai_family() == socket.AF_INET


if os.environ.get("OTSO_IPV4", "").strip() not in ("", "0"):
    prefer_ipv4()


def check_ipv6(url: str) -> None:
    """Once per host: switch to IPv4 if the host has an IPv6 address that does not answer."""
    parts = urlsplit(url)
    host = parts.hostname
    port = parts.port or (443 if parts.scheme == "https" else 80)
    if not host or host in _PROBED or ipv4_only():
        return
    _PROBED.add(host)
    try:
        v6 = socket.getaddrinfo(host, port, socket.AF_INET6, socket.SOCK_STREAM)
        socket.getaddrinfo(host, port, socket.AF_INET, socket.SOCK_STREAM)
    except OSError:
        return  # no IPv6 address (or no IPv4 one): nothing to choose between
    try:
        with socket.socket(socket.AF_INET6, socket.SOCK_STREAM) as s:
            s.settimeout(IPV6_PROBE_TIMEOUT)
            s.connect(v6[0][4])
    except OSError:
        prefer_ipv4()


def _timeout(timeout):
    if timeout is None or isinstance(timeout, tuple):
        return timeout
    return (min(CONNECT_TIMEOUT, timeout), timeout)


def get(url: str, *, timeout=30, **kwargs) -> requests.Response:
    """requests.get with the IPv6 check and a short connect timeout."""
    check_ipv6(url)
    return requests.get(url, timeout=_timeout(timeout), **kwargs)


def head(url: str, *, timeout=10, **kwargs) -> requests.Response:
    """requests.head with the IPv6 check and a short connect timeout."""
    check_ipv6(url)
    return requests.head(url, timeout=_timeout(timeout), **kwargs)
