import uuid

import pytest

from config import get_test_dir, get_ws_client


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

 
