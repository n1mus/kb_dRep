import uuid

import pytest

from config import get_test_dir, get_ws_client
from data import *


@pytest.fixture
def test_dir():
    return get_test_dir()

@pytest.fixture(scope='session')
def ws():
    ws_name = "kb_dRep_" + str(uuid.uuid4())
    ws_client = get_ws_client()
    ws_id = ws_client.create_workspace({'workspace': ws_name})[0]                      
    yield {                                                                           
        'workspace_id': ws_id,                                                               
        'workspace_name': ws_name,                                                           
    } 
    ws_client.delete_workspace({'workspace': ws_name})

@pytest.fixture
def kb_clients():
    kb_clients = {'dfu': get_mock_dfu(), 'mgu': get_mock_mgu(), 'au': get_mock_au()}
    return kb_clients


