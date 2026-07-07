import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "backend"))
import uvicorn
if __name__ == "__main__": uvicorn.run("qtlift.api:app", host="127.0.0.1", port=8765, reload=False, app_dir=str(Path(__file__).parent/"backend"))

