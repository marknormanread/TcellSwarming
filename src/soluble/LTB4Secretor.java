package soluble;

import java.util.Vector;

import core.SecretionRecord;

public interface LTB4Secretor extends Secretor
{
	public boolean isSecretingLTB4();
	
	public double getLTB4SecretionRate();
	
	public Vector<SecretionRecord> getSecretionHistory();
}
