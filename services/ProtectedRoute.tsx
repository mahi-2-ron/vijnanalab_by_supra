import { Navigate } from 'react-router-dom';
import { useAuth } from '../services/AuthContext';

type Props = {
  children: JSX.Element;
  allowedRoles?: string[];
};

const ProtectedRoute = ({ children, allowedRoles }: Props) => {
  const { user, role, loading } = useAuth();

  if (loading) {
    return (
      <div className="h-screen flex items-center justify-center text-white">
        Loading...
      </div>
    );
  }

  if (!user) {
    return <Navigate to="/login" replace />;
  }

  if (allowedRoles && !allowedRoles.includes(role || '')) {
    return <Navigate to="/home" replace />;
  }

  return children;
};

export default ProtectedRoute;
